# timestep_routine.py
# Module for timestep handling

import numpy as np
import time
import pylab as pl

from PyLCM.parameters import *
from PyLCM.micro_particle import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from PyLCM.collision import *
from PyLCM.animation import *
from PyLCM.widget import *

from Post_process.analysis import *
from Post_process.print_plot import *

def timesteps_function(n_particles_widget, P_widget, RH_widget, T_widget, w_widget, nt_widget, dt_widget, rm_spec, ascending_mode_widget, mode_displaytype_widget, z_widget, max_z_widget, Condensation_widget, Collision_widget, mode_aero_init_widget, gridwidget, kohler_activation_radius, switch_kappa_koehler, switch_sedi_removal): 

    
    # Function call of the complete model initialization (model_init) (aerosol initialization included)
    P_parcel, T_parcel, q_parcel, z_parcel, w_parcel, N_aero, mu_aero, sigma_aero, nt, dt, \
    max_z, do_condensation, do_collision, ascending_mode, time_half_wave_parcel, S_lst, display_mode, \
    qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, \
    z_parcel_array, particles_list, spectra_arr, con_ts, act_ts, evp_ts, dea_ts, acc_ts, aut_ts, precip_ts, particles_array, rc_liq_avg_array, rc_liq_std_array, particles_c_array,n_particles = model_init(dt_widget, nt_widget, Condensation_widget, Collision_widget, \
                                n_particles_widget, T_widget, P_widget, RH_widget, w_widget, z_widget, \
                                max_z_widget, mode_aero_init_widget, gridwidget, \
                                ascending_mode_widget, mode_displaytype_widget,switch_kappa_koehler)  

################################
    # Timestep routine
################################

    # Create array for the drop radii evolution
    

    if display_mode == 'graphics':
        # Initialization of animation
        figure_item = animation_init(dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array)

    for t in range(nt):
        time = (t+1)*dt
        # Parcel ascending
        z_parcel, T_parcel, rho_parcel, V_parcel, air_mass_parcel = ascend_parcel(z_parcel, T_parcel,P_parcel, w_parcel, dt, time,max_z, time_half_wave_parcel, ascending_mode)
        
        # Condensational Growth
        dq_liq = 0.0
        if do_condensation:
            particles_list, T_parcel, q_parcel, S_lst, con_ts[t+1], act_ts[t+1], evp_ts[t+1], dea_ts[t+1] = drop_condensation(particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass_parcel, S_lst, rho_aero,kohler_activation_radius, con_ts[t+1], act_ts[t+1], evp_ts[t+1], dea_ts[t+1],switch_kappa_koehler)
            
            # Convert mass output to per mass per sec.
            con_ts[t+1]  = 1e3 * con_ts[t+1] / air_mass_parcel / dt
            act_ts[t+1]  = 1e3 * act_ts[t+1] / air_mass_parcel / dt
            evp_ts[t+1]  = 1e3 * evp_ts[t+1] / air_mass_parcel / dt
            dea_ts[t+1]  = 1e3 * dea_ts[t+1] / air_mass_parcel / dt
            
        # Collisional Growth
        if do_collision:
            particles_list, acc_ts[t+1], aut_ts[t+1],precip_ts[t+1] = collection(dt, particles_list,rho_parcel, rho_liq, P_parcel, T_parcel, acc_ts[t+1], aut_ts[t+1],precip_ts[t+1], switch_sedi_removal, z_parcel, max_z, w_parcel)
            
            # Convert mass output to per mass per sec.
            acc_ts[t+1]  = 1e3 * acc_ts[t+1] / air_mass_parcel / dt
            aut_ts[t+1]  = 1e3 * aut_ts[t+1] / air_mass_parcel / dt
        # Analysis
        spectra_arr[t+1],qa_ts[t+1], qc_ts[t+1],qr_ts[t+1], na_ts[t+1], nc_ts[t+1], nr_ts[t+1], particles_array[t+1], rc_liq_avg_array[t+1], rc_liq_std_array[t+1], particles_c_array[t+1] = ts_analysis(particles_list,air_mass_parcel,rm_spec, n_bins, n_particles)
        RH_parcel = (q_parcel * P_parcel / (q_parcel + r_a / rv)) / esatw( T_parcel ) 
        
        # Saving values of T_parcel, RH_parcel, q_parcel, z_parcel for every timestep (needed for plots)
        T_parcel_array[t+1]  = T_parcel
        RH_parcel_array[t+1] = RH_parcel
        q_parcel_array[t+1]  = q_parcel
        z_parcel_array[t+1]  = z_parcel
                
        time_array = np.arange(nt+1)*dt

        # Display of variables during runtime
        if display_mode == 'text_fast':
            # Prints text output at every second
            if (time%1) ==0:
                print_output(t,dt, z_parcel, T_parcel, q_parcel, RH_parcel, qc_ts[t+1], qr_ts[t+1], na_ts[t+1], nc_ts[t+1], nr_ts[t+1])
        elif display_mode == 'graphics':
            # Displays and continuously updates plots during runtime using plotly
            # Figure output is updated every 5 seconds
            if (time%5) == 0:
                animation_call(figure_item, time_array, t, dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array)
            

    return nt, dt, time_array, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, qa_ts,qc_ts,qr_ts, na_ts,nc_ts,nr_ts, spectra_arr, con_ts, act_ts, evp_ts, dea_ts, acc_ts, aut_ts, particles_array, rc_liq_avg_array, rc_liq_std_array, particles_c_array
    
