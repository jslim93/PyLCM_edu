from PyLCM.parameters import *
from PyLCM.micro_particle import *
import numpy as np
from PyLCM.parcel import *
from PyLCM.condensation import *
from Post_process.print_plot import *

from scipy.stats import lognorm

def model_init(dt_widget, nt_widget, Condensation_widget, Collision_widget, n_particles_widget, T_widget, P_widget, RH_widget, w_widget, z_widget, max_z_widget, mode_aero_init_widget, gridwidget, ascending_mode_widget, mode_displaytype_widget):
    # Data of the model steering parameters is being read out from the widgets.
    # Values necessary for the model initialization will be returned.
    
    # Model steering parameters
    dt = dt_widget.value # Default: 0.5 s
    nt = nt_widget.value # Default: 3600

    do_condensation = Condensation_widget.value  # Default: True
    do_collision    = Collision_widget.value  # Default: False

    n_particles = n_particles_widget.value

    # Parcel initial parameters
    T_parcel   = T_widget.value
    P_parcel   = P_widget.value
    RH_parcel  = RH_widget.value
    w_parcel   = w_widget.value
    z_parcel   = z_widget.value

    # RH to q conversion
    q_parcel    = RH_parcel * esatw( T_parcel ) / ( P_parcel - RH_parcel * esatw( T_parcel ) ) * r_a / rv
        
    max_z = max_z_widget.value
    
    # Aerosol initialization
    mode_aero_init = mode_aero_init_widget.value  # Possible values are: "weighting_factor", 'random'

    # Read in the variables given above taking into account unit factors.
    N_aero     = [gridwidget[1, 0].value*1.0E6, gridwidget[1, 1].value*1.0E6, gridwidget[1, 2].value*1.0E6, gridwidget[1, 3].value*1.0E6]
    mu_aero    = [gridwidget[2, 0].value*1.0E-6, gridwidget[2, 1].value*1.0E-6, gridwidget[2, 2].value*1.0E-6, gridwidget[2, 3].value*1.0E-6]
    sigma_aero = [gridwidget[3, 0].value, gridwidget[3, 1].value, gridwidget[3, 2].value, gridwidget[3, 3].value]
    
    # Truncate the array before taking the log if one of the N_aero_i is 0, which means that this mode will no longer be used.
    N_aero_array = np.array(N_aero) # First: Convert into np.array
    zeroindices  = np.where(N_aero_array==0) # Get the number of the mode which is empty
    zeroindices  = zeroindices[0]       # Conversion for better usage

    # Conversion of the other indices
    mu_aero_array = np.array(mu_aero)
    sigma_aero_array = np.array(sigma_aero)

    # Now delete the respective item in each array (N, mu, sigma)
    if len(zeroindices) > 0:
        # delete
        N_aero_array     = np.delete(N_aero_array, zeroindices)
        mu_aero_array    = np.delete(mu_aero_array, zeroindices)
        sigma_aero_array = np.delete(sigma_aero_array, zeroindices)

    # Now perform the log of mu and sigma arrays
    mu_aero_array = np.log(mu_aero_array)
    sigma_aero_array = np.log(sigma_aero_array)

    # Renaming of the arrrays
    N_aero = N_aero_array
    mu_aero = mu_aero_array
    sigma_aero = sigma_aero_array
    
    # Further initialization
    dz=0
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)
    
    # Aerosol initialization
    T_parcel, q_parcel, particles_list = aero_init(mode_aero_init, n_particles, P_parcel, T_parcel,q_parcel, N_aero, mu_aero, sigma_aero, rho_aero)
    
    # Parcel routine
    # Initalize spectrum output
    spectra_arr = np.zeros((nt+1,len(rm_spec)))
    # Initialization of arrays for time series output
    qa_ts,qc_ts,qr_ts = np.zeros(nt+1),np.zeros(nt+1),np.zeros(nt+1)
    na_ts,nc_ts,nr_ts = np.zeros(nt+1),np.zeros(nt+1),np.zeros(nt+1)
    con_ts, act_ts, evp_ts, dea_ts = np.zeros(nt+1),np.zeros(nt+1),np.zeros(nt+1),np.zeros(nt+1)
    acc_ts, aut_ts = np.zeros(nt+1),np.zeros(nt+1)
    spectra_arr[0],qa_ts[0], qc_ts[0],qr_ts[0], na_ts[0], nc_ts[0], nr_ts[0] = ts_analysis(particles_list,air_mass_parcel,rm_spec, n_bins)
    
    
    # Initialization of arrays for T_parcel, RH_parcel, q_parcel and z_parcel. 
    # They will later be filled with values for each time step.
    T_parcel_array  = np.zeros(nt+1)
    RH_parcel_array = np.zeros(nt+1)
    q_parcel_array  = np.zeros(nt+1)
    z_parcel_array  = np.zeros(nt+1)


    # Inserting the initialization values at the 0th position of the arrays.
    T_parcel_array[0]  = T_parcel
    RH_parcel_array[0] = (q_parcel * P_parcel / (q_parcel + r_a / rv)) / esatw( T_parcel ) 
    q_parcel_array[0]  = q_parcel
    z_parcel_array[0]  = z_parcel
    
    # Read in the ascending mode selected in the widget
    ascending_mode=ascending_mode_widget.value
    # Settings for the 'sine' and the 'in_cloud_oscillation' modes: time half wavelength of the parcel (s)
    time_half_wave_parcel = 600.0  # This value can be adapted by the user

    S_lst = 0.0
    
    # Read in selected display mode (options are: 'text_fast' or 'graphics')
    display_mode = mode_displaytype_widget.value
    
    return P_parcel, T_parcel, q_parcel, z_parcel, w_parcel, N_aero, mu_aero, sigma_aero, nt, dt, max_z, do_condensation, do_collision, ascending_mode, time_half_wave_parcel, S_lst, display_mode, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, particles_list, spectra_arr, con_ts, act_ts, evp_ts, dea_ts, acc_ts, aut_ts



def aero_init(mode_aero_init, n_ptcl, P_parcel, T_parcel,q_parcel, N_aero, mu_aero,sigma_aero,rho_aero):
    # Aerosol initialization
    
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)
    
    particles_list = []
    # Computation of saturated water vapour pressure (e_s) and water vapour pressure of the parcel (e_a)
    e_s     = 611.2 * np.exp( 17.62 * ( T_parcel - 273.15 ) / ( T_parcel - 29.65 ) )
    e_a     = q_parcel * P_parcel / ( q_parcel + r_a / rv )
    # Computation of supersaturation
    S_adia  = ( e_a - e_s ) / e_s

    dql_liq = np.sum([p.M for p in particles_list])

    min_mass_aero = 1.0E-200 
    
    # Computation of aerosol inititial radius
    radius = np.logspace(np.log10(1.0E-9), np.log10(1.0E-6), n_ptcl)
    dlogr   = ( np.log(2.0E-6) - np.log(1.0E-9) ) / n_ptcl
        
    # Initialize particle (aerosol particles)
    for i in range(n_ptcl):
        particle = particles(i)
        
        # RNG (random number generator) method to generate distribution where all particles represent the same number of droplets
        if mode_aero_init == "random":
            mode_count = len(N_aero)

            n_particles_mode_array = np.zeros(mode_count) # Initialization of an array for n_particles of the 4 modes
            
            # Assigning the number of droplets to the 4 modes according to specified n_particles and N_aero[i]
            # (The result is now no longer integer!)

            for j in range(mode_count):
                n_particles_mode_array[j] = N_aero[j] / np.sum(N_aero) * n_ptcl
            
            # Truncate particle modes to integer and sum them up
            n_particles_mode_int = n_particles_mode_array.astype(int)
            n_particles_intsum = np.sum(n_particles_mode_int) # Array with integer values of particle numbers
            
            # Sum of particles of different modes (float)
            n_particles_floatsum = np.sum(n_particles_mode_array)
            
            # Compare the integer and the float sum
            n_difference = n_particles_floatsum - n_particles_intsum
            
            # Add the missing droplets (or remove the droplets which are too many) to the first mode
            n_particles_mode_int[0] = int(n_particles_mode_int[0] + n_difference)
            
            # Generate log-normal distribution for the modes
            for k in range(mode_count):
                if k == 0:
                    temp_arr = np.random.lognormal(mu_aero[0], sigma_aero[0], n_particles_mode_int[0])
                else:
                    temp_arr  =  np.concatenate((temp_arr, np.random.lognormal(mu_aero[k], sigma_aero[k], n_particles_mode_int[k])))
            aero_r_seed = temp_arr
                       
            particle.A = air_mass_parcel * np.sum(N_aero)/n_ptcl
            particle.Ns = aero_r_seed[i]**3 * 4./3. * np.pi * rho_aero * particle.A
        
        # Bin-like method to generate distribution where each particle represents different number of droplets
        elif mode_aero_init == "weighting_factor":
            
            # Calculate the PDF (probability density function) of the overlapping log-normal distributions
            pdf_sum = np.zeros_like(radius)
            for N, mu, sigma in zip(N_aero, mu_aero, sigma_aero):
                pdf = lognormal_pdf(radius, mu,sigma)
                pdf_sum += N * pdf
            # Define the range of values to evaluate the PDF
            particle.A = air_mass_parcel * pdf_sum[i] * dlogr * radius[i]
            particle.Ns = radius[i]**3 * 4./3. * np.pi * rho_aero * particle.A

        particles_list.append(particle)

        if particle.Ns > min_mass_aero:
            r_aero = (particle.Ns/ ( particle.A * 4.0 / 3.0 * pi * rho_aero ) )**(1.0/3.0)

            particle.M = max(r_aero, r_equi(S_adia,T_parcel,r_aero, rho_aero))**3 * particle.A * 4.0 / 3.0 * pi * rho_liq
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