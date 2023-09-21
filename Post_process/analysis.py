import numpy as np
from matplotlib import pyplot as plt
import time
import pylab as pl
import pandas as pd
from IPython import display
from PyLCM.parameters import *
from PyLCM.micro_particle import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from PyLCM.collision import *

 
def ts_analysis(particles_list,air_mass_parcel,log_edges, nbins):
    # Timesteps analysis: Performs calculations of q_x mixing ratios and n_x number densities and the spectra
    nbins = nbins - 1 # Number of bins are 1 smaller than number of edges.
    
    spec = np.zeros(nbins)
    # Calculate the total mass of particles with r_liq larger than activation_radius_ts
    qc_mass = 0.0
    qr_mass = 0.0
    qa_mass = 0.0
    
    NC = 0.0
    NR = 0.0
    NA = 0.0

    for particle in particles_list:
        r_liq = (particle.M / (particle.A * 4.0 / 3.0 * np.pi * rho_liq))**(1/3.0)
        if r_liq > activation_radius_ts:
            if r_liq < seperation_radius_ts:
                qc_mass += particle.M
                NC += particle.A
            else:
                qr_mass += particle.M
                NR += particle.A
        else:
            qa_mass += particle.M 
            NA += particle.A
            
        spec = get_spec(nbins,spec,log_edges,r_liq,particle.A,air_mass_parcel)

    qc = qc_mass / air_mass_parcel *1e3
    qr = qr_mass / air_mass_parcel *1e3
    qa = qa_mass / air_mass_parcel *1e3
    
    NA = NA / air_mass_parcel /1e6
    NC = NC / air_mass_parcel /1e6
    NR = NR / air_mass_parcel /1e6
    
    return(spec,qa, qc,qr, NA, NC, NR)

def get_spec(nbins,spectra_arr,log_edges,r_liq,weight_factor,air_mass_parcel):
    # Computes array of the spectra
    bin_idx = np.searchsorted(log_edges, r_liq, side='right') - 1
    if 0 <= bin_idx < nbins:
        spectra_arr[bin_idx] += weight_factor / air_mass_parcel / (rr_spec[bin_idx] - rl_spec[bin_idx])*rm_spec[bin_idx]

    return spectra_arr

def save_model_output_variables(time_array, RH_parcel_array, q_parcel_array, T_parcel_array, z_parcel_array, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, filename='testoutput_model.csv'):
    # Saves the output arrays to a csv file in the subfolder 'Output'
    # Optional filename can be given
    
    output_variables_array = np.stack((time_array, RH_parcel_array, q_parcel_array*1000, T_parcel_array, z_parcel_array, qa_ts*1000, qc_ts*1000, qr_ts*1000, na_ts/1e6, nc_ts/1e6, nr_ts/1e6), axis=-1)
    
    # Conversion to pandas format
    output_variables_dataframe = pd.DataFrame(output_variables_array)
    output_variables_dataframe.columns=['time', 'RH_parcel', 'q_parcel', 'T_parcel', 'z_parcel', 'qa_ts', 'qc_ts', 'qr_ts', 'na_ts', 'nc_ts', 'nr_ts']
    
    # Save to csv
    output_variables_dataframe.to_csv('Output/'+filename)
    print('Output data written to: Output/'+filename)
    
def save_model_output_dsd(spectra_arr, rm_spec, rl_spec, rr_spec, nt, filename='dsd_array_output.csv'):
    # Saves the output of the droplet size distributions to a csv-file. The filename can be adjusted manually.
    # The first 3 columns contain radii of the bin edges and bin mean
    # Remaining columns: each column representing the DSD of one timestep (timestep-number as defined in column heading)
    
    # Create list of timesteps for row names
    timesteplist = list(range(nt+1))
    firstnames = ['rm_spec [µm]', 'rl_spec [µm]', 'rr_spec [µm]']
    # Combine the lists => list of all row names
    rowlist = firstnames + timesteplist
    
    # Attatch the columns of rm, rl, rr to the spectra array
    dsd_array = np.column_stack((rm_spec*1e6, rl_spec*1e6, rr_spec*1e6, spectra_arr.T/1e6))
    # Division /1e6 in spectra_arr as done for the DSD plots, to match the unit cm⁻3
    dsd_array = dsd_array.T
    # Convert to pandas DataFrame and assigning rowlist as index column (row names)
    dsd_dataframe = pd.DataFrame(dsd_array, index=rowlist)

    # Save to csv
    dsd_dataframe.to_csv('Output/'+filename)
    print('Output data of droplet size distribution written to: Output/'+filename)   