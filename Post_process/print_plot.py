import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors
import copy
from IPython.display import clear_output
from matplotlib import cm

from PyLCM.parameters import *
from PyLCM.micro_init import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from PyLCM.collision import *
from Post_process.analysis import *


def spec_plot(ax, spectra_arr, nt,dt, rm_spec):
    # Clear the previous plot
    clear_output(wait=True)

    # Set the y-scale to logarithmic
    ax.set_yscale("log")

    # Create the contour plot with a logarithmic color scale
    contour = ax.contourf(np.arange(nt+1)*dt, 1e6*rm_spec, spectra_arr.T, norm=matplotlib.colors.LogNorm())

    # Create a logarithmic colorbar
    cbar = plt.colorbar(contour, ax=ax, format="%.0e")

    # Add labels and title
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radius [µm]")
    ax.set_title("DSD time evolution")
    ax.set_ylim([1e-2,1e4])    
    
    
def print_output(t,dt, z_parcel, T_parcel, q_parcel, rh, qc, qr, na, nc, nr):
    # Clear previous output
    clear_output(wait=True)

    # Print the initial variable names
    print("value: {:<8}  {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}".format(
        "Time (s)", "z (m)", "T (K)", "qv (g/kg)", "RH (%)", "QC (g/kg)", "QR (g/kg)", "NA (/mg)", "NC (/mg)", "NR (/mg)"))

    # Print the updated output
    print("after: {:<8.1f}  {:<8.2f} {:<8.2f} {:<9.2f} {:<8.3f}  {:<8.3f}  {:<8.3f}  {:<8.2f}  {:<8.2f}  {:<8.2f}".format(
        (t+1) * dt, z_parcel, T_parcel, 1e3 * q_parcel, 100* rh, 1e3 * qc, 1e3 * qr, na / 1e6, nc / 1e6, nr / 1e6))
    
def subplot_array_function(plot_mode, dt, nt, rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, spectra_arr, increment_widget, con_ts, act_ts, evp_ts, dea_ts, acc_ts, aut_ts):
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
    # plot 1: DSD size distribution
    spec_plot(axs[1,0],spectra_arr/1e6, nt,dt,rm_spec)
    
    # plot 2: particle densities
    # calculate number of drawn spectra (nt_spec) via nt of the model and user given increment (via increment_widget)
    line_increment = increment_widget.value
    nt_spec = nt / line_increment
    nt_spec = int(nt_spec)
    cmap = cm.get_cmap('jet')
    norm = plt.Normalize(0, nt_spec - 1)
        
    # spectras in user given timesteps (defined via line_increment) will be displayed
    for i in range(nt_spec):
        spectra_arr_nan = copy.deepcopy(spectra_arr) 
        # deepcopy needed so that values <= 0 are masked out only in the copy for the plot (and remain for the data output)
        spectra_arr_nan[np.where(spectra_arr_nan<=0)] = np.nan
        axs[1,1].plot(rm_spec*1e6, spectra_arr_nan[i*line_increment]/1e6, color=cmap(norm(i)))
        axs[1,1].set_yscale("log")
        axs[1,1].set_xscale("log")
        axs[1,1].set_xlabel('radius [µm]')
        axs[1,1].set_ylabel('particle densities dN/dlog(R) [mg$^{-1}$]')
        #axs[1,1].set_ylim(1)
        
    # add colorbar for particle densities plot
    # select the plot axis
    plotaxis = axs[1,1]
    # get the colormap
    cmap_spectra = plt.cm.get_cmap('jet')
    # normalize the colormap to the max. timestep, take into account the timestep interval dt [s]
    norm2 = plt.Normalize(vmin=0, vmax=np.max((nt_spec*line_increment-line_increment)*dt))
    # produce mappable object
    scalarmap = plt.cm.ScalarMappable(norm=norm2, cmap=cmap_spectra)
    # add colorbar
    fig.colorbar(scalarmap, ax=plotaxis, orientation='vertical', label='Time [t]')

    # plot 3: plot for condensation and evaporation
    if plot_mode=='time-series':
        axs[1,2].plot(time_array, con_ts, label = "Condensation", color='darkblue')
        axs[1,2].plot(time_array, act_ts, label = "Activation", color='limegreen', linestyle=':')
        axs[1,2].plot(time_array, evp_ts, label = "Evaporaton", color='brown')
        axs[1,2].plot(time_array, dea_ts, label = "Deactivation", color='black', linestyle='--')
        axs[1,2].set_xlabel("Time [s]")
        axs[1,2].set_ylabel("Conversion Rates [kg$^{-1}$s$^{-1}$]")
        axs[1,2].legend()
        
    elif plot_mode=='vertical profile':
        axs[1,2].plot(con_ts, z_parcel_array, label = "Condensation", color='darkblue')
        axs[1,2].plot(act_ts, z_parcel_array, label = "Activation", color='limegreen', linestyle=':')
        axs[1,2].plot(evp_ts, z_parcel_array, label = "Evaporaton", color='brown')
        axs[1,2].plot(dea_ts, z_parcel_array, label = "Deactivation", color='black', linestyle='--')
        axs[1,2].set_xlabel("Conversion Rates [kg$^{-1}$s$^{-1}$]")
        axs[1,2].set_ylabel("Height $z$ [m]")
        axs[1,2].legend()
        
    
    # plot 4: autoconversion and accretion
    # compute moving average of both arrays and the time_array, z_parcel_array using built in Pandas routine
    # set the moving window length
    window_length = 120
    # moving average and first window_length-1 steps deleted (which are nan)
    aut_ts_rolling = pd.Series(aut_ts).rolling(window=window_length).mean().iloc[window_length-1:].values
    acc_ts_rolling = pd.Series(acc_ts).rolling(window=window_length).mean().iloc[window_length-1:].values
    time_rolling = pd.Series(time_array).rolling(window=window_length).mean().iloc[window_length-1:].values
    z_parcel_rolling = pd.Series(z_parcel_array).rolling(window=window_length).mean().iloc[window_length-1:].values
    
    if plot_mode=='time-series':
        axs[1,3].plot(time_rolling, aut_ts_rolling, label = "Autoconversion", color='purple')
        axs[1,3].plot(time_rolling, acc_ts_rolling, label = "Accretion", color='darkorange', linestyle='--')
        axs[1,3].set_xlabel("Time [s]")
        axs[1,3].set_ylabel("Conversion Rates [kg$^{-1}$s$^{-1}$] \n rolling mean, window size: "+str(window_length))
        axs[1,3].legend()
        
    elif plot_mode=='vertical profile':
        axs[1,3].plot(aut_ts_rolling, z_parcel_rolling, label = "Autoconversion", color='purple')
        axs[1,3].plot(acc_ts_rolling, z_parcel_rolling, label = "Accretion", color='darkorange', linestyle='--')
        axs[1,3].set_xlabel("Conversion Rates [kg$^{-1}$s$^{-1}$] \n rolling mean, window size: "+str(window_length))
        axs[1,3].set_ylabel("Height $z$ [m]")
        axs[1,3].legend()
    
    
    
    fig.tight_layout()
    fig.show()
