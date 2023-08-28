import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors
import copy
from IPython.display import clear_output
from matplotlib import cm
from parameters import *
from micro import *
from aero_init import *
from parcel import *
from condensation import *
from collision import *
from analysis import *

# module for animated plots using plotly
# version 1 for 4 plots (testing)
# help about plotly: https://plotly.com/python/creating-and-updating-figures/
# https://plotly.com/python/line-charts/
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import time
from IPython.display import clear_output


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
    
def subplot_array_function(plot_mode, dt, nt, rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, spectra_arr):
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
    spec_plot(axs[1,0],spectra_arr/1e6, nt,dt,rm_spec)
    nt_spec = 180
    cmap = cm.get_cmap('jet')
    norm = plt.Normalize(0, nt_spec - 1)
    # particle densities
    for i in range(nt_spec):
        spectra_arr_nan = copy.deepcopy(spectra_arr) 
        # deepcopy needed so that values <= 0 are masked out only in the copy for the plot (and remain for the data output)
        spectra_arr_nan[np.where(spectra_arr_nan<=0)] = np.nan
        axs[1,1].plot(rm_spec*1e6, spectra_arr_nan[i*20]/1e6, color=cmap(norm(i)))
        axs[1,1].set_yscale("log")
        axs[1,1].set_xscale("log")
        axs[1,1].set_xlabel('radius [µm]')
        axs[1,1].set_ylabel('particle densities dN/dlog(R) [mg$^{-1}$]')
        #axs[1,1].set_ylim(1)

    fig.tight_layout()
    fig.show()

####################################################################################################################################################################################
#Animation Routine
####################################################################################################################################################################################


def animation_init(dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array):
    # sets all initial parameters for the animation
    time_array = np.arange(nt+1)*dt
    
    # Create subplot layout
    fig = make_subplots(rows=3, cols=2)

    # Add initial traces (lines) to subplots
    fig.add_trace(go.Scatter(x=time_array, y=RH_parcel_array*100, mode='lines', line_color='lightblue', name='RH [%]', showlegend=True), row=1, col=1)
    fig.add_trace(go.Scatter(x=time_array, y=q_parcel_array*1e3, mode='lines', line_color='lightgreen', name='q_v [g/kg]'), row=1, col=2)
    fig.add_trace(go.Scatter(x=time_array, y=z_parcel_array, mode='lines', line_color='black', name='z [m]'), row=2, col=1)
    fig.add_trace(go.Scatter(x=time_array, y=T_parcel_array, mode='lines', line_color='red', name='T [K]'), row=2, col=2)
    # Plot for the mixing ratios (with 3 lines: aerosol, cloud, rain)
    fig.add_trace(go.Scatter(x=time_array, y=qa_ts*1e3, mode='lines', line_color='blue', name='q_a [g/kg] Aerosol'), row=3, col=1)
    fig.add_trace(go.Scatter(x=time_array, y=qc_ts*1e3, mode='lines', line_color='orange', name='q_c [g/kg] Cloud'), row=3, col=1)
    fig.add_trace(go.Scatter(x=time_array, y=qr_ts*1e3, mode='lines', line_color='green', name='q_r [g/kg] Rain'), row=3, col=1)
    # Plot for the number concentrations (not yet set to log axis)
    fig.add_trace(go.Scatter(x=time_array, y=na_ts/1e6, mode='lines', line_color='blue', name='n_a [mg⁻1] Aerosol'), row=3, col=2)
    fig.add_trace(go.Scatter(x=time_array, y=nc_ts/1e6, mode='lines', line_color='orange', name='n_c [mg⁻1] Cloud'), row=3, col=2)
    fig.add_trace(go.Scatter(x=time_array, y=nr_ts/1e6, mode='lines', line_color='green', name='n_r [mg⁻1] Rain'), row=3, col=2)
    
    
    # some layout commands
    fig.update_layout(height=850, width=1250, title='Time series', yaxis_title='Relative Humidity RH [-]', yaxis_range=[RH_parcel_array[0]*100,104], yaxis2_title='q_v [g/kg]', yaxis3_title='Height z [m]', xaxis5_title='Time [s]', xaxis6_title='Time [s]', yaxis4_title='Temperature T [K]')
    # axis title for the last 2 subplots
    fig.update_layout(yaxis5_title='Liquid Water Mixing ratios q_x [g/kg]', yaxis6_title='Number concentrations n_x [mg⁻1]')
    # axis title font size changed to lower value
    fig.update_layout(yaxis_title_font_size=11, yaxis1_title_font_size=11, yaxis2_title_font_size=11, yaxis3_title_font_size=11, yaxis4_title_font_size=11, yaxis5_title_font_size=11, yaxis6_title_font_size=11)
    
    
    # fix x axis range
    fig.update_layout(xaxis_range=[0, (nt+1)*dt], xaxis2_range=[0, (nt+1)*dt], xaxis3_range=[0, (nt+1)*dt], xaxis4_range=[0, (nt+1)*dt], xaxis5_range=[0, (nt+1)*dt], xaxis6_range=[0, (nt+1)*dt])
    # change y axis of last subplot to log scale
    fig.update_layout(yaxis6_type="log")
    
    
    return fig
    
def animation_call(fig, time_array, t, dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array):
    # to be included inside the loop
    
    # Clear previous figures
    clear_output(wait=True)
    
    # Update traces (lines) in all subplots
    fig.data[0].x = time_array[0:t]
    fig.data[0].y = RH_parcel_array*100
    fig.data[1].x = time_array[0:t]
    fig.data[1].y = q_parcel_array*1e3
    fig.data[2].x = time_array[0:t]
    fig.data[2].y = z_parcel_array
    fig.data[3].x = time_array[0:t]
    fig.data[3].y = T_parcel_array
    # Update traces for the mixing ratios plot
    fig.data[4].x = time_array[0:t]
    fig.data[4].y = qa_ts*1e3
    fig.data[5].x = time_array[0:t]
    fig.data[5].y = qc_ts*1e3
    fig.data[6].x = time_array[0:t]
    fig.data[6].y = qr_ts*1e3
    # Update traces for the number concentrations plot
    fig.data[7].x = time_array[0:t]
    fig.data[7].y = na_ts/1e6
    fig.data[8].x = time_array[0:t]
    fig.data[8].y = nc_ts/1e6
    fig.data[9].x = time_array[0:t]
    fig.data[9].y = nr_ts/1e6
    
    
    # fig.update_layout for the title
    #static_title = 'Time series: Current model output values: '
    print_output_str = 'Time (s): '+str((t+1)*dt)+', z (m): '+str(np.round(z_parcel_array[t], decimals=2))+', RH (%): '+str(np.round(RH_parcel_array[t]*100, decimals=3))+', T (K): '+str(np.round(T_parcel_array[t], decimals=2))
    current_title = print_output_str
    fig.update_layout(title=current_title)
    
    # Display updated figure
    fig.show()

    # Pause for a while
    # time.sleep(0.0001)
