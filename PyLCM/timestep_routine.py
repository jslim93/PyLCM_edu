# timestep_routine.py
# module for timestep handling

import numpy as np
import time
import pylab as pl

from PyLCM.parameters import *
from PyLCM.micro import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from PyLCM.collision import *
from PyLCM.animation import *
from PyLCM.widget import *

from Post_process.analysis import *
from Post_process.print_plot import *

def timesteps_function(mode_aero_init, n_particles, P_parcel, T_parcel,q_parcel, z_parcel, w_parcel, N_aero, mu_aero,\
                                                                       sigma_aero, rho_aero, molecular_weight_aero, nt, dt, rm_spec, \
                                                                       ascending_mode_widget, mode_displaytype_widget, max_z, \
                                                                       do_condensation, do_collision, kohler_activation_radius, act_crit_r):
    dz=0
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)
    #Aerosol init
    T_parcel, q_parcel, particles_list = aero_init(mode_aero_init, n_particles, P_parcel, T_parcel,q_parcel, N_aero, mu_aero, sigma_aero, rho_aero, molecular_weight_aero)
    #parcel routine
    #initalize spectrum output
    spectra_arr = np.zeros((nt+1,len(rm_spec)))
    # init of array for time series output
    qa_ts,qc_ts,qr_ts = np.zeros(nt+1),np.zeros(nt+1),np.zeros(nt+1)
    na_ts,nc_ts,nr_ts = np.zeros(nt+1),np.zeros(nt+1),np.zeros(nt+1)
    con_ts, act_ts, evp_ts, dea_ts = 0.0, 0.0, 0.0, 0.0
    spectra_arr[0],qa_ts[0], qc_ts[0],qr_ts[0], na_ts[0], nc_ts[0], nr_ts[0] = qc_qr_analysis(particles_list,air_mass_parcel,rm_spec, n_bins)
    
    # init of array for T_parcel, RH_parcel, q_parcel and z_parcel values for each timestep
    T_parcel_array  = np.zeros(nt+1)
    RH_parcel_array = np.zeros(nt+1)
    q_parcel_array  = np.zeros(nt+1)
    z_parcel_array  = np.zeros(nt+1)


    # inserting the init. values to the 0th position of the arrays
    T_parcel_array[0]  = T_parcel
    RH_parcel_array[0] = (q_parcel * P_parcel / (q_parcel + r_a / rv)) / esatw( T_parcel ) 
    q_parcel_array[0]  = q_parcel
    z_parcel_array[0]  = z_parcel
    
    ascending_mode=ascending_mode_widget.value
    time_half_wave_parcel = 600.0  # maybe change to widget or variable input later

    S_lst = 0.0
    
    # read in display mode
    display_mode = mode_displaytype_widget.value
    
    if display_mode == 'graphics':
        # initialization of animation
        figure_item = animation_init(dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array)

    for t in range(nt):
        time = (t+1)*dt
        #Parcel ascending
        z_parcel, T_parcel, rho_parcel, V_parcel, air_mass_parcel = ascend_parcel(z_parcel, T_parcel,P_parcel, w_parcel, dt, time, time_half_wave_parcel, ascending_mode, max_z)
        
        #Condensational Growth
        dq_liq = 0.0
        if do_condensation:
                
            particles_list, T_parcel, q_parcel, S_lst, con_ts, act_ts, evp_ts, dea_ts = drop_condensation(particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass_parcel, S_lst, rho_aero,kohler_activation_radius, act_crit_r,con_ts, act_ts, evp_ts, dea_ts)
        #Collisional Growth
        if do_collision:
            particles_list = collection(dt, particles_list,rho_parcel, rho_liq, P_parcel, T_parcel)
            
        #Analysis
        spectra_arr[t+1],qa_ts[t+1], qc_ts[t+1],qr_ts[t+1], na_ts[t+1], nc_ts[t+1], nr_ts[t+1] = qc_qr_analysis(particles_list,air_mass_parcel,rm_spec, n_bins)
        RH_parcel = (q_parcel * P_parcel / (q_parcel + r_a / rv)) / esatw( T_parcel ) 
        #convert mass out put to per mass every 1 sec.
        if (time%1) ==0:
            con_ts /= air_mass_parcel
            act_ts /= air_mass_parcel
            evp_ts /= air_mass_parcel
            dea_ts /= air_mass_parcel
            
            con_ts, act_ts, evp_ts, dea_ts = 0.0, 0.0, 0.0, 0.0
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
            

    return time_array, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, qa_ts,qc_ts,qr_ts, na_ts,nc_ts,nr_ts, spectra_arr
    
