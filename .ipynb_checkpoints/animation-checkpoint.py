# module for animated plots using plotly
# version 1 for 4 plots (testing)
# help about plotly: https://plotly.com/python/creating-and-updating-figures/
# https://plotly.com/python/line-charts/
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import time
from IPython.display import clear_output
import numpy as np


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
