import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors
from IPython.display import clear_output
from parameters import *
from micro import *
from aero_init import *
from parcel import *
from condensation import *
from collision import *
from analysis import *

def spec_plot(ax, spectra_arr, nt,dt, log_edges):
    # Clear the previous plot
    clear_output(wait=True)

    # Set the y-scale to logarithmic
    ax.set_yscale("log")

    # Create the contour plot with a logarithmic color scale
    contour = ax.contourf(np.arange(nt+1)*dt, log_edges[1:], spectra_arr.T, norm=matplotlib.colors.LogNorm())

    # Create a logarithmic colorbar
    cbar = plt.colorbar(contour, ax=ax, format="%.0e")

    # Add labels and title
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radius [µm]")
    ax.set_title("DSD time evolution")
    ax.set_ylim([1e-8,2e-3])
    
    
    
    
def print_output(t,dt, z_parcel, T_parcel, q_parcel, rh, qc, qr, na, nc, nr):
    # Clear previous output
    clear_output(wait=True)

    # Print the initial variable names
    print("value: {:<8}  {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}".format(
        "Time (s)", "z (m)", "T (K)", "qv (g/kg)", "RH (%)", "QC (g/kg)", "QR (g/kg)", "NA (/mg)", "NC (/mg)", "NR (/mg)"))

    # Print the updated output
    print("after: {:<8.1f}  {:<8.2f} {:<8.2f} {:<9.2f} {:<8.3f}  {:<8.3f}  {:<8.3f}  {:<8.2f}  {:<8.2f}  {:<8.2f}".format(
        (t+1) * dt, z_parcel, T_parcel, 1e3 * q_parcel, 100* rh, 1e3 * qc, 1e3 * qr, na / 1e6, nc / 1e6, nr / 1e6))
    
def subplot_array_function(plot_mode, dt, nt, log_edges, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, spectra_arr):
    # t not needed?
    # initialisation of subplot layout
    fig, axs = plt.subplots(2, 4, sharex=False, sharey=False, figsize=(18,8))
    
    time_array = np.arange(nt+1)*dt


    # 1st row
    # left: q_v and RH
    if plot_mode=='time-series':
        ax2 = axs[0,0].twinx() # secondary y axis
        axs[0,0].plot(time_array, RH_parcel_array*100, label = "RH [%]")
        axs[0,0].axhline(y=100, color='gray', linestyle='-', linewidth=0.5) # horizontal line at RH=100%
        ax2.plot(time_array, q_parcel_array*1e3, label = "$q_v$ [g/kg]", color='green')
        axs[0,0].set_xlabel("Time [s]")
        axs[0,0].set_ylabel('relative humidity RH [%]')
        ax2.set_ylabel('$q_v$ [g/kg]')
        #axs[0,0].set_ylim([70,110])
        axs[0,0].legend(loc='lower right')
        ax2.legend(loc='lower center')
    elif plot_mode=='vertical profile':
        ax2 = axs[0,0].twiny() # secondary x axis, shared y axis (height z)
        axs[0,0].plot(RH_parcel_array*100, z_parcel_array, label = "RH [%]")
        axs[0,0].axvline(x=100, color='gray', linestyle='-', linewidth=0.5) # vertical line at RH=100%
        axs[0,0].set_xlabel("relative humidity RH [%]")
        axs[0,0].set_ylabel('Height $z$ [m]')
        ax2.plot(q_parcel_array*1e3, z_parcel_array, label = "$q_v$ [g/kg]", color='green')
        ax2.set_xlabel('$q_v$ [g/kg]')
        axs[0,0].legend(loc='lower right')
        ax2.legend(loc='lower center')
    
    
    # middle: T
    if plot_mode=='time-series':
        ax2_2 = axs[0,1].twinx() # secondary y axis for this plot
        axs[0,1].plot(time_array, T_parcel_array, label = "$T$ [K]", color='red')
        ax2_2.plot(time_array, z_parcel_array, label = "$z$ [m]", color='black')
        axs[0,1].set_xlabel("Time [s]")
        axs[0,1].set_ylabel('temperature $T$ [K]')
        ax2_2.set_ylabel('height $z$ [m]')
        axs[0,1].legend(loc='lower right')
        ax2_2.legend(loc='lower center')
    elif plot_mode=='vertical profile':
        axs[0,1].plot(T_parcel_array, z_parcel_array, label = "T [K]", color='red')
        axs[0,1].set_xlabel('Temperature $T$ [K]')
        axs[0,1].set_ylabel('Height $z$ [m]')
        axs[0,1].grid()

  
    # mixing ratios
    if plot_mode=='time-series':
        axs[0,2].plot(time_array,qa_ts*1e3, label = "Aerosol")
        axs[0,2].plot(time_array,qc_ts*1e3, label = "Cloud")
        axs[0,2].plot(time_array,qr_ts*1e3, label = "Rain")
        axs[0,2].set_xlabel("Time [s]")
        axs[0,2].set_ylabel('mixing ratios $q_x$ [g/kg]')
        axs[0,2].legend()
        axs[0,2].grid()
    elif plot_mode=='vertical profile':
        axs[0,2].plot(qa_ts*1e3, z_parcel_array, label = "Aerosol")
        axs[0,2].plot(qc_ts*1e3, z_parcel_array, label = "Cloud")
        axs[0,2].plot(qr_ts*1e3, z_parcel_array, label = "Rain")
        axs[0,2].set_ylabel("Height $z$ [m]")
        axs[0,2].set_xlabel('mixing ratios $q_x$ [g/kg]')
        axs[0,2].legend()
        axs[0,2].grid()

    # right: number concentrations
    if plot_mode=='time-series':
        axs[0,3].plot(time_array,na_ts/1e6, label = "Aerosol")
        axs[0,3].plot(time_array,nc_ts/1e6, label = "Cloud")
        axs[0,3].plot(time_array,nr_ts/1e6, label = "Rain")
        axs[0,3].set_xlabel("Time [s]")
        axs[0,3].set_ylabel('number concentration of particles $n_x$ [mg$^{-1}$]')
        axs[0,3].set_yscale('log')
        axs[0,3].legend()
        axs[0,3].grid()
    elif plot_mode=='vertical profile':
        axs[0,3].plot(na_ts/1e6, z_parcel_array, label = "Aerosol")
        axs[0,3].plot(nc_ts/1e6, z_parcel_array, label = "Cloud")
        axs[0,3].plot(nr_ts/1e6, z_parcel_array, label = "Rain")
        axs[0,3].set_ylabel("Height $z$ [m]")
        axs[0,3].set_xlabel('number concentration of particles $n_x$ [mg$^{-1}$]')
        axs[0,3].set_xscale('log')
        axs[0,3].grid()


    # 2nd row
    # DSD size distribution
    spec_plot(axs[1,0],spectra_arr/1e6, nt,dt,log_edges)

    # particle densities
    for i in range(18):
        axs[1,1].plot(log_edges[1:]*1e6, spectra_arr[i*100]/1e6)
        axs[1,1].set_yscale("log")
        axs[1,1].set_xscale("log")
        axs[1,1].set_xlabel('radius [µm]')
        axs[1,1].set_ylabel('particle densities N [cm$^{-3}$]')
        axs[1,1].set_ylim([1e-1,100])

    fig.tight_layout()