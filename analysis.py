import numpy as np
from matplotlib import pyplot as plt
import time
import pylab as pl
import pandas as pd
from IPython import display
from parameters import *
from micro import *
from aero_init import *
from parcel import *
from condensation import *
from collision import *

def qc_qr_analysis(particles_list,air_mass_parcel,log_edges, nbins):
    nbins = nbins - 1 # number of bins are 1 smaller than number of edges.
    
    spec = np.zeros(nbins)
     
    activation_radius_ts = 1.0E-6  # Activation radius in meters, example value
    seperation_radius_ts = 25.0E-6
    # Calculate the total mass of particles with r_liq larger than activation_radius_ts
    qc_mass = 0.0
    qr_mass = 0.0
    NC = 0.0
    NR = 0.0
    qa = 0.0
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
            qa += particle.M / air_mass_parcel
            NA += particle.A / air_mass_parcel
            
        spec = get_spec(nbins,spec,log_edges,r_liq,particle.A,air_mass_parcel)

    qc = qc_mass  / air_mass_parcel
    qr = qr_mass  / air_mass_parcel
 
    NC = NC  / air_mass_parcel
    NR = NR  / air_mass_parcel
    
    return(spec,qa, qc,qr, NA, NC, NR)

def get_spec(nbins,spectra_arr,log_edges,r_liq,weight_factor,air_mass_parcel):
    from parameters import  rr_spec
    from parameters import  rl_spec
    from parameters import  rm_spec
    
    bin_idx = np.searchsorted(log_edges, r_liq, side='right') - 1
    if 0 <= bin_idx < nbins:
        spectra_arr[bin_idx] += weight_factor / air_mass_parcel / (rr_spec[bin_idx] - rl_spec[bin_idx])*rm_spec[bin_idx]
        

    return spectra_arr

def save_model_output_variables(time_array, RH_parcel_array, q_parcel_array, T_parcel_array, z_parcel_array, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, filename='testoutput_model.csv'):
    # saves the output arrays to a csv file in the subfolder 'output'
    # optional filename can be given
    
    output_variables_array = np.stack((time_array, RH_parcel_array, q_parcel_array*1000, T_parcel_array, z_parcel_array, qa_ts*1000, qc_ts*1000, qr_ts*1000, na_ts/1e6, nc_ts/1e6, nr_ts/1e6), axis=-1)
    
    # conversion to pandas
    output_variables_dataframe = pd.DataFrame(output_variables_array)
    output_variables_dataframe.columns=['time', 'RH_parcel', 'q_parcel', 'T_parcel', 'z_parcel', 'qa_ts', 'qc_ts', 'qr_ts', 'na_ts', 'nc_ts', 'nr_ts']
    
    # save to csv
    output_variables_dataframe.to_csv('output/'+filename)
    print('Output data written to: output/'+filename)
    
def save_model_output_dsd(spectra_arr, rm_spec, rl_spec, rr_spec, nt, filename='dsd_array_output.csv'):
    # saves the output of the droplet size distributions to a csv-file, filename can be adjusted manually
    # first 3 columns: radii of the bin edges and bin mean
    # remaining columns: each column representing the DSD of one timestep (timestep-number like in column heading)
    
    # create list of timesteps for row names
    timesteplist = list(range(nt+1))
    firstnames = ['rm_spec [µm]', 'rl_spec [µm]', 'rr_spec [µm]']
    # combine the lists => list of all row names
    rowlist = firstnames + timesteplist
    
    # attatch the columns of rm, rl, rr to the spectra array
    dsd_array = np.column_stack((rm_spec*1e6, rl_spec*1e6, rr_spec*1e6, spectra_arr.T/1e6))
    # new: division /1e6 in spectra_arr as done for the DSD plots, to match the unit cm⁻3
    dsd_array = dsd_array.T
    # convert to pandas DataFrame and assigning rowlist as index column (row names)
    dsd_dataframe = pd.DataFrame(dsd_array, index=rowlist)

    # save to csv
    dsd_dataframe.to_csv('output/'+filename)
    print('Output data of droplet size distribution written to: output/'+filename)   
    
    
