# timestep_routine.py
# module for timestep handling

import numpy as np
import time
import pylab as pl
from IPython import display
from parameters import *
from micro import *
from aero_init import *
from parcel import *
from condensation import *
from collision import *
from analysis import *
from print_plot import *
from animation import *
from widget import *

def timesteps_function(n_particles_widget, P_widget, RH_widget, T_widget, w_widget, nt_widget, dt_widget, rm_spec, ascending_mode_widget, mode_displaytype_widget, max_z_widget, Condensation_widget, Collision_widget, mode_aero_init_widget, gridwidget):
    
      
    # call of the complete model initialization (model_init) (aerosol initialization included)
    P_parcel, T_parcel, q_parcel, z_parcel, w_parcel, N_aero, mu_aero, sigma_aero, nt, dt, \
    max_z, do_condensation, do_collision, ascending_mode, time_half_wave_parcel, S_lst, display_mode, \
    qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, \
    z_parcel_array, particles_list, spectra_arr = model_init(dt_widget, nt_widget, Condensation_widget, Collision_widget, \
                                n_particles_widget, T_widget, P_widget, RH_widget, w_widget, \
                                max_z_widget, mode_aero_init_widget, gridwidget, \
                                ascending_mode_widget, mode_displaytype_widget)
    
        
    ########
    # timestep routine
    ########
    
    if display_mode == 'graphics':
        # initialization of animation
        figure_item = animation_init(dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array)


    for t in range(nt):
        time = (t+1)*dt
        #Parcel ascending
        #if z_parcel < max_z: 
        z_parcel, T_parcel, rho_parcel, V_parcel, air_mass_parcel = ascend_parcel(z_parcel, T_parcel,P_parcel, w_parcel, dt, time, time_half_wave_parcel, ascending_mode, max_z=max_z)
        
        #Condensational Growth
        dq_liq = 0.0
        if do_condensation:
            particles_list, T_parcel, q_parcel, S_lst = drop_condensation(particles_list, T_parcel, q_parcel, P_parcel, dt, air_mass_parcel,S_lst, rho_aero, molecular_weight_aero)

        #Collisional Growth
        if do_collision:
            particles_list = collection(dt, particles_list,rho_parcel, rho_liq, P_parcel, T_parcel)
            
        #Analysis
        spectra_arr[t+1],qa_ts[t+1], qc_ts[t+1],qr_ts[t+1], na_ts[t+1], nc_ts[t+1], nr_ts[t+1] = qc_qr_analysis(particles_list,air_mass_parcel,rm_spec, n_bins)
        RH_parcel = (q_parcel * P_parcel / (q_parcel + r_a / rv)) / esatw( T_parcel ) 
        
        # saving of T_parcel, RH_parcel, q_parcel, z_parcel for every timestep (needed for plots)
        T_parcel_array[t+1]  = T_parcel
        RH_parcel_array[t+1] = RH_parcel
        q_parcel_array[t+1]  = q_parcel
        z_parcel_array[t+1]  = z_parcel
        
        time_array = np.arange(nt+1)*dt

        if display_mode == 'text_fast':
            #Visulaization at every second
            if (time%1) ==0:
                print_output(t,dt, z_parcel, T_parcel, q_parcel, RH_parcel, qc_ts[t+1], qr_ts[t+1], na_ts[t+1], nc_ts[t+1], nr_ts[t+1])
        elif display_mode == 'graphics':
            # function to draw and update the plotly figure, every 5 seconds (idea for time saving)
            if (time%5) == 0:
                animation_call(figure_item, time_array, t, dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array)
            

    return time_array, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, qa_ts,qc_ts,qr_ts, na_ts,nc_ts,nr_ts, spectra_arr, dt, nt
    
