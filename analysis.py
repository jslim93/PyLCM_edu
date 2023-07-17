import numpy as np
from matplotlib import pyplot as plt
import time
import pylab as pl
from IPython import display
from parameters import *
from micro import *
from chem import *
from aero_init import *
from parcel import *
from condensation import *
from collision import *

def qc_qr_analysis(particles_list,air_mass_parcel,log_edges):
    nbins = 100
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
    
    bin_idx = np.searchsorted(log_edges, r_liq, side='right') - 1
    if 0 <= bin_idx < nbins:
        spectra_arr[bin_idx] += weight_factor / air_mass_parcel

    return spectra_arr